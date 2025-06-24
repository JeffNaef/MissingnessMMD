import torch
import torch.nn as nn
from statistics import mean
import matplotlib.pyplot as plt
from types import SimpleNamespace

# Make sure these local files are in the same directory or Python path
import loss 
from model import Generator, Discriminator
from data import NoiseGenerator, generate_contaminated_data, NoEndingDataLoaderIter
from torch.utils.data import TensorDataset
from train import train_one_round
from utils import set_seed, initialize_d_optimizer, initialize_g_optimizer, coord_median, plot_visualization

def run(num_epoch=150, num_iter=-1, p=100, s=-1,
        sparse_estimation=True, eps=0.2, train_size=50000,
        coord_median_as_origin=True, contamination="gauss_5",
        loss_name="JSLoss", kappa=None, l1_constrain_type="scale",
        real_batch_size=500, fake_batch_size=500, debug=False,
        simultaneous=True, num_step_d=1, num_step_g=1,
        d_optimizer="adam", g_optimizer="sgd",
        d_sgd_lr=0.02, d_sgd_momentum=0.9, sgd_weight_decay=0,
        d_adam_lr=0.0002, d_adam_b1=0.5, d_adam_b2=0.999,
        adam_weight_decay=0,
        d_adagrad_lr=0.01, d_adagrad_lr_decay=0,
        d_adagrad_initial_accumulator_value=0.0,
        adagrad_weight_decay=0,
        g_sgd_lr=0.02, g_sgd_momentum=0.0,
        g_adam_lr=0.0002, g_adam_b1=0.5, g_adam_b2=0.999,
        real_grad_penalty=None, fake_grad_penalty=None,
        seed=0):
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    
    # Use SimpleNamespace to allow attribute access (args.p instead of args['p'])
    args = SimpleNamespace(
        num_epoch=num_epoch, num_iter=num_iter, p=p, s=s,
        sparse_estimation=sparse_estimation, eps=eps, train_size=train_size,
        coord_median_as_origin=coord_median_as_origin, contamination=contamination,
        loss=loss_name, kappa=kappa, l1_constrain_type=l1_constrain_type,
        real_batch_size=real_batch_size, fake_batch_size=fake_batch_size,
        debug=debug, simultaneous=simultaneous, num_step_d=num_step_d,
        num_step_g=num_step_g, d_optimizer=d_optimizer, g_optimizer=g_optimizer,
        d_sgd_lr=d_sgd_lr, d_sgd_momentum=d_sgd_momentum, sgd_weight_decay=sgd_weight_decay,
        d_adam_lr=d_adam_lr, d_adam_b1=d_adam_b1, d_adam_b2=d_adam_b2,
        adam_weight_decay=adam_weight_decay, d_adagrad_lr=d_adagrad_lr,
        d_adagrad_lr_decay=d_adagrad_lr_decay,
        d_adagrad_initial_accumulator_value=d_adagrad_initial_accumulator_value,
        adagrad_weight_decay=adagrad_weight_decay, g_sgd_lr=g_sgd_lr,
        g_sgd_momentum=g_sgd_momentum, g_adam_lr=g_adam_lr, g_adam_b1=g_adam_b1,
        g_adam_b2=g_adam_b2, real_grad_penalty=real_grad_penalty,
        fake_grad_penalty=fake_grad_penalty, seed=seed
    )


    print(args)

    if args.s == -1:
        args.s = None

    assert not (args.num_epoch == -1 and args.num_iter == -1)
    # assert args.real_batch_size <= args.train_size
    if args.debug and args.p != 2:
        raise ValueError(args.debug, args.p)

    assert isinstance(args.kappa, (int, float, dict, type(None)))
    if isinstance(args.kappa, dict):
        init_kappa = args.kappa[0]
    else:
        init_kappa = args.kappa

    set_seed(args.seed)

    if args.s is None:
        theta = torch.zeros(args.p).to(device)
    else:
        theta = torch.zeros(args.p).to(device)
        theta[0:args.s] = 1.

    data, theta = generate_contaminated_data(
        args.eps, args.train_size,
        theta=theta,
        type_cont=args.contamination,
        coord_median_as_origin=args.coord_median_as_origin)
        
        
        
    data = data.to(device)
    theta = theta.to(device)

    data_loader = torch.utils.data.DataLoader(
        TensorDataset(data),
        batch_size=args.real_batch_size, shuffle=True, num_workers=0)

    lst_activation = [nn.Sigmoid()]
    lst_num_hidden = [20]

    loss_obj = getattr(loss, args.loss)()
    noise_generator = NoiseGenerator().to(device)
    generator = Generator(
        p=args.p,
        initializer=coord_median(data_loader.dataset.tensors[0]),
        # initializer = torch.ones(args.p) * 2
    ).to(device)
    discriminator = Discriminator(
        input_dim=args.p,
        lst_num_hidden=lst_num_hidden,
        lst_activation=lst_activation,
        kappa=init_kappa,
        l1_constrain_type=args.l1_constrain_type,
    ).to(device)

    d_optim = initialize_d_optimizer(discriminator.parameters(), args)
    g_optim = initialize_g_optimizer(generator.parameters(), args)

    print("dist {:.4f}".format(torch.norm(generator.eta - theta).item()))

    data_loader_iter = NoEndingDataLoaderIter(data_loader)


    # g_scheduler = torch.optim.lr_scheduler.StepLR(
    #     g_optim, step_size=1, gamma=0.98)
    # d_scheduler = torch.optim.lr_scheduler.StepLR(
    #     d_optim, step_size=1, gamma=0.98)

    epoch = 0
    idx_iter = 0

    lst_eta = [generator.get_numpy_eta()]
    while True:
        idx_iter += 1

        if isinstance(args.kappa, dict):
            if epoch in args.kappa.keys():
                discriminator.kappa = args.kappa[epoch]
                print("Set kappa to {}".format(discriminator.kappa))
                del args.kappa[epoch]

        # XXX: note that training does not stop exactly at the end of the epoch

        lst_d_loss, lst_g_loss = train_one_round(
            loss_obj, discriminator, generator, d_optim, g_optim,
            data_loader_iter, noise_generator,
            fake_batch_size=args.fake_batch_size,
            device=None,
            real_grad_penalty=args.real_grad_penalty,
            fake_grad_penalty=args.fake_grad_penalty,
            num_step_d=args.num_step_d, num_step_g=args.num_step_g,
            simultaneous=args.simultaneous,
            s=args.s,
            sparse_estimation=args.sparse_estimation,
        )

        if data_loader_iter.epoch > epoch:
            lst_eta.append(generator.get_numpy_eta())
            print(
                "epoch {:6d},".format(epoch),
                "dist {:.4f},".format(
                    torch.norm(generator.eta - theta).item()),
                "d_loss {:.4f},".format(mean(lst_d_loss)),
                "g_loss {:.4f},".format(mean(lst_g_loss)),
            )
            epoch = data_loader_iter.epoch

            if args.num_epoch != -1 and \
                    data_loader_iter.epoch >= args.num_epoch:
                break

            if args.debug:
                fig = plt.figure()
                fig.set_size_inches((10, 8))
                plot_visualization(
                    discriminator, generator, data_loader, theta,
                    device=None)
                title = 'Epoch ' + str(epoch)
                plt.title(title, fontsize=15)
                fig.savefig('./Figure/' + title + '.png')
                plt.close()

                print(generator.get_numpy_eta())

        if args.num_iter != -1 and idx_iter > args.num_iter:
            break

    # Return the results instead of saving to a file
    print("Training complete.")
    return (theta.cpu().numpy(), lst_eta)
